/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifdef WIN32
#include <direct.h>
#endif

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>

#ifdef  _FINE_GRAIN_MPI
#include <mpi.h>
#endif



#ifdef _USE_PTHREADS
#include <pthread.h>

#endif

#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#include <xmmintrin.h>
/*
  special bug fix, enforces denormalized numbers to be flushed to zero,
  without this program is a tiny bit faster though.
  #include <emmintrin.h> 
  #define MM_DAZ_MASK    0x0040
  #define MM_DAZ_ON    0x0040
  #define MM_DAZ_OFF    0x0000
*/
#endif

#include "axml.h"
#include "globalVariables.h"


#define _PORTABLE_PTHREADS


/***************** UTILITY FUNCTIONS **************************/


void myBinFwrite(const void *ptr, size_t size, size_t nmemb)
{ 
  size_t  
    bytes_written = fwrite(ptr, size, nmemb, byteFile);

  assert(bytes_written == nmemb);
}





void *malloc_aligned(size_t size) 
{
  void 
    *ptr = (void *)NULL;
 
  int 
    res;
  

#if defined (__APPLE__)
  /* 
     presumably malloc on MACs always returns 
     a 16-byte aligned pointer
  */

  ptr = malloc(size);
  
  if(ptr == (void*)NULL) 
   assert(0);
  
#ifdef __AVX
  assert(0);
#endif


#else
  res = posix_memalign( &ptr, BYTE_ALIGNMENT, size );

  if(res != 0) 
    assert(0);
#endif 
   
  return ptr;
}








void printBothOpen(const char* format, ... )
{
  FILE *f = myfopen(infoFileName, "ab");

  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);

  fclose(f);
}

void printBothOpenMPI(const char* format, ... )
{
#ifdef _WAYNE_MPI
  if(processID == 0)
#endif
    {
      FILE *f = myfopen(infoFileName, "ab");

      va_list args;
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
      
      va_start(args, format);
      vprintf(format, args );
      va_end(args);
      
      fclose(f);
    }
}


boolean getSmoothFreqs(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}

const unsigned int *getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}


int getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].states;
}

int getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].undetermined;
}



char getInverseMeaning(int dataType, unsigned char state)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return  pLengths[dataType].inverseMeaning[state];
}

partitionLengths *getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
    states    = p->states,
    tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  pLength.leftLength = pLength.rightLength = states * states;
  pLength.eignLength = states;
  pLength.evLength   = states * states;
  pLength.eiLength   = states * states;
  pLength.substRatesLength = (states * states - states) / 2;
  pLength.frequenciesLength = states;
  pLength.tipVectorLength   = tipLength * states;
  pLength.symmetryVectorLength = (states * states - states) / 2;
  pLength.frequencyGroupingLength = states;
  pLength.nonGTR = FALSE;

  return (&pLengths[dataType]); 
}







double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}

int gettimeSrand(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 24*60*60*localtm.tm_yday + 60*60*localtm.tm_hour + 60*localtm.tm_min  + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec;
#endif
}

double randum (long  *seed)
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}

static int filexists(char *filename)
{
  FILE *fp;
  int res;
  fp = fopen(filename,"rb");

  if(fp)
    {
      res = 1;
      fclose(fp);
    }
  else
    res = 0;

  return res;
}


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("\n Error: the file %s you want to open for reading does not exist, exiting ...\n\n", path);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }
  else
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("\n Error: the file %s you want to open for writing or appending can not be opened [mode: %s], exiting ...\n\n",
		   path, mode);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }


}





/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/










/***********************reading and initializing input ******************/

static void getnums (rawdata *rdta)
{
  if (fscanf(INFILE, "%d %d", & rdta->numsp, & rdta->sites) != 2)
    {
      if(processID == 0)
	printf("\n Error: problem reading number of species and sites\n\n");
      errorExit(-1);
    }

  if (rdta->numsp < 4)
    {
      if(processID == 0)
	printf("\n Error: too few species\n\n");
      errorExit(-1);
    }

  if (rdta->sites < 1)
    {
      if(processID == 0)
	printf("\n Error: too few sites\n\n");
      errorExit(-1);
    }

  return;
}





boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}


static void uppercase (int *chptr)
{
  int  ch;

  ch = *chptr;
  if ((ch >= 'a' && ch <= 'i') || (ch >= 'j' && ch <= 'r')
      || (ch >= 's' && ch <= 'z'))
    *chptr = ch + 'A' - 'a';
}




static void getyspace (rawdata *rdta)
{
  size_t size = 4 * ((size_t)(rdta->sites / 4 + 1));
  int    i;
  unsigned char *y0;

  rdta->y = (unsigned char **) malloc((rdta->numsp + 1) * sizeof(unsigned char *));
  assert(rdta->y);   

  y0 = (unsigned char *) malloc(((size_t)(rdta->numsp + 1)) * size * sizeof(unsigned char));

  /*
    printf("Raw alignment data Assigning %Zu bytes\n", ((size_t)(rdta->numsp + 1)) * size * sizeof(unsigned char));

  */

  assert(y0);   

  rdta->y0 = y0;

  for (i = 0; i <= rdta->numsp; i++)
    {
      rdta->y[i] = y0;
      y0 += size;
    }

  return;
}


static unsigned int KISS32(void)
{
  static unsigned int 
    x = 123456789, 
    y = 362436069,
    z = 21288629,
    w = 14921776,
    c = 0;

  unsigned int t;

  x += 545925293;
  y ^= (y<<13); 
  y ^= (y>>17); 
  y ^= (y<<5);
  t = z + w + c; 
  z = w; 
  c = (t>>31); 
  w = t & 2147483647;

  return (x+y+w);
}

static boolean setupTree (tree *tr, analdef *adef)
{
  nodeptr  p0, p, q;
  int
    i,
    j,   
    tips,
    inter; 

  if(!adef->readTaxaOnly)
    {
      /*tr->bigCutoff = FALSE;*/

      tr->patternPosition = (int*)NULL;
      tr->columnPosition = (int*)NULL;

      /*tr->maxCategories = MAX(4, adef->categories);*/

      /*tr->partitionContributions = (double *)malloc(sizeof(double) * tr->NumberOfModels);

      for(i = 0; i < tr->NumberOfModels; i++)
	tr->partitionContributions[i] = -1.0;

      tr->perPartitionLH = (double *)malloc(sizeof(double) * tr->NumberOfModels);
      

      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  tr->perPartitionLH[i] = 0.0;	 
	}

      if(adef->grouping)
	tr->grouped = TRUE;
      else
	tr->grouped = FALSE;

      if(adef->constraint)
	tr->constrained = TRUE;
      else
	tr->constrained = FALSE;

	tr->treeID = 0;*/
    }

  tips  = tr->mxtips;
  inter = tr->mxtips - 1;

  if(!adef->readTaxaOnly)
    {
      tr->yVector      = (unsigned char **)  malloc((tr->mxtips + 1) * sizeof(unsigned char *));

      /*      tr->fracchanges  = (double *)malloc(tr->NumberOfModels * sizeof(double));
	      tr->likelihoods  = (double *)malloc(adef->multipleRuns * sizeof(double));*/
    }

  /*tr->numberOfTrees = -1;

 

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)calloc(tr->treeStringLength, sizeof(char)); 
  tr->tree0 = (char*)calloc(tr->treeStringLength, sizeof(char));
  tr->tree1 = (char*)calloc(tr->treeStringLength, sizeof(char));*/


  /*TODO, must that be so long ?*/

  if(!adef->readTaxaOnly)
    {
            
      /*tr->td[0].count = 0;
      tr->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * tr->mxtips);
      tr->td[0].executeModel = (boolean *)malloc(sizeof(boolean) * tr->NumberOfModels);
      tr->td[0].parameterValues = (double *)malloc(sizeof(double) * tr->NumberOfModels);
       
      for(i = 0; i < tr->NumberOfModels; i++)
	tr->fracchanges[i] = -1.0;
      tr->fracchange = -1.0;

      tr->constraintVector = (int *)malloc((2 * tr->mxtips) * sizeof(int));*/

      tr->nameList = (char **)malloc(sizeof(char *) * (tips + 1));
    }

  if (!(p0 = (nodeptr) malloc((tips + 3*inter) * sizeof(node))))
    {
      printf("\n Error: unable to obtain sufficient tree memory\n\n");
      return  FALSE;
    }
  
  /*  tr->nodeBaseAddress = p0;


  if (!(tr->nodep = (nodeptr *) malloc((2*tr->mxtips) * sizeof(nodeptr))))
    {
      printf("\n Error: unable to obtain sufficient tree memory, too\n\n");
      return  FALSE;
      }*/

  //tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  

  tr->vLength = 0;

  tr->h = (hashtable*)NULL;


  return TRUE;
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

static boolean getdata(analdef *adef, rawdata *rdta, tree *tr)
{
  int   
    i, 
    j, 
    basesread, 
    basesnew, 
    ch, my_i, meaning,
    len,
    meaningAA[256], 
    meaningDNA[256], 
    meaningBINARY[256],
    meaningGeneric32[256],
    meaningGeneric64[256];
  
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

  for (i = 0; i < 256; i++)
    {      
      meaningAA[i]          = -1;
      meaningDNA[i]         = -1;
      meaningBINARY[i]      = -1;
      meaningGeneric32[i]   = -1;
      meaningGeneric64[i]   = -1;
    }

  /* generic 32 data */

  for(i = 0; i < 32; i++)
    meaningGeneric32[genericChars32[i]] = i;
  meaningGeneric32['-'] = getUndetermined(GENERIC_32);
  meaningGeneric32['?'] = getUndetermined(GENERIC_32);

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

  meaningAA['X'] = 
    meaningAA['?'] = 
    meaningAA['*'] = 
    meaningAA['-'] = 
    getUndetermined(AA_DATA);

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

  meaningDNA['N'] = 
    meaningDNA['O'] = 
    meaningDNA['X'] = 
    meaningDNA['-'] = 
    meaningDNA['?'] = 
    getUndetermined(DNA_DATA);

  /* BINARY DATA */

  meaningBINARY['0'] = 1;
  meaningBINARY['1'] = 2;
  
  meaningBINARY['-'] = 
    meaningBINARY['?'] = 
    getUndetermined(BINARY_DATA);


  /*******************************************************************/

  basesread = basesnew = 0;

  allread = FALSE;
  firstpass = TRUE;
  ch = ' ';

  while (! allread)
    {
      for (i = 1; i <= tr->mxtips; i++)
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
		      if(processID == 0)
			{
			  printf("Taxon Name to long at taxon %d, adapt constant nmlngth in\n", i);
			  printf("axml.h, current setting %d\n", nmlngth);
			}
		      errorExit(-1);
		    }
		}
	      while(ch !=  ' ' && ch != '\n' && ch != '\t' && ch != '\r');

	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
		ch = getc(INFILE);
	      
	      ungetc(ch, INFILE);

	      buffer[my_i] = '\0';
	      len = strlen(buffer) + 1;
	      checkTaxonName(buffer, len);
	      tr->nameList[i] = (char *)malloc(sizeof(char) * len);
	      strcpy(tr->nameList[i], buffer);
	    }

	  j = basesread;

	  while ((j < rdta->sites) && ((ch = getc(INFILE)) != EOF) && (ch != '\n') && (ch != '\r'))
	    {
	      uppercase(& ch);

	      assert(tr->dataVector[j + 1] != -1);

	      switch(tr->dataVector[j + 1])
		{
		case BINARY_DATA:
		  meaning = meaningBINARY[ch];
		  break;
		case DNA_DATA:
		case SECONDARY_DATA:
		case SECONDARY_DATA_6:
		case SECONDARY_DATA_7:
		  /*
		     still dealing with DNA/RNA here, hence just act if as they where DNA characters
		     corresponding column merging for sec struct models will take place later
		  */
		  meaning = meaningDNA[ch];
		  break;
		case AA_DATA:
		  meaning = meaningAA[ch];
		  break;
		case GENERIC_32:
		  meaning = meaningGeneric32[ch];
		  break;
		case GENERIC_64:
		  meaning = meaningGeneric64[ch];
		  break;
		default:
		  assert(0);
		}

	      if (meaning != -1)
		{
		  j++;
		  rdta->y[i][j] = ch;		 
		}
	      else
		{
		  if(!whitechar(ch))
		    {
		      printf("\n Error: bad base (%c) at site %d of sequence %d\n\n",
			     ch, j + 1, i);
		      return FALSE;
		    }
		}
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
			   j - basesread, basesnew - basesread, i, tr->nameList[i]);
		    return  FALSE;
		  }
	    }
	  while (ch != '\n' && ch != EOF && ch != '\r') ch = getc(INFILE);  /* flush line *//* PC-LINEBREAK*/
	}

      firstpass = FALSE;
      basesread = basesnew;
      allread = (basesread >= rdta->sites);
    }

  for(j = 1; j <= tr->mxtips; j++)
    for(i = 1; i <= rdta->sites; i++)
      {
	assert(tr->dataVector[i] != -1);

	switch(tr->dataVector[i])
	  {
	  case BINARY_DATA:
	    meaning = meaningBINARY[rdta->y[j][i]];
	    if(meaning == getUndetermined(BINARY_DATA))
	      gaps++;
	    break;

	  case SECONDARY_DATA:
	  case SECONDARY_DATA_6:
	  case SECONDARY_DATA_7:
	    assert(tr->secondaryStructurePairs[i - 1] != -1);
	    assert(i - 1 == tr->secondaryStructurePairs[tr->secondaryStructurePairs[i - 1]]);
	    /*
	       don't worry too much about undetermined column count here for sec-struct, just count
	       DNA/RNA gaps here and worry about the rest later-on, falling through to DNA again :-)
	    */
	  case DNA_DATA:
	    meaning = meaningDNA[rdta->y[j][i]];
	    if(meaning == getUndetermined(DNA_DATA))
	      gaps++;
	    break;

	  case AA_DATA:
	    meaning = meaningAA[rdta->y[j][i]];
	    if(meaning == getUndetermined(AA_DATA))
	      gaps++;
	    break;

	  case GENERIC_32:
	    meaning = meaningGeneric32[rdta->y[j][i]];
	    if(meaning == getUndetermined(GENERIC_32))
	      gaps++;
	    break;

	  case GENERIC_64:
	    meaning = meaningGeneric64[rdta->y[j][i]];
	    if(meaning == getUndetermined(GENERIC_64))
	      gaps++;
	    break;
	  default:
	    assert(0);
	  }

	total++;
	rdta->y[j][i] = meaning;
      }

  adef->gapyness = (double)gaps / (double)total;
    
  /*myBinFwrite(&(adef->gapyness), sizeof(double), 1);*/

  printf("gappyness: %f\n", adef->gapyness);
  
  /*for(i = 1; i <= tr->mxtips; i++)
    {
      int 
	len = strlen(tr->nameList[i]) + 1;
      
      myBinFwrite(&len, sizeof(int), 1);
      myBinFwrite(tr->nameList[i], sizeof(char), len);
      
      printf("%d %s\n", len, tr->nameList[i]);
      }     */
  
  return  TRUE;
}



static void inputweights (rawdata *rdta)
{
  int i, w, fres;
  FILE *weightFile;
  int *wv = (int *)malloc(sizeof(int) *  rdta->sites);

  weightFile = myfopen(weightFileName, "rb");

  i = 0;

  while((fres = fscanf(weightFile,"%d", &w)) != EOF)
    {
      if(!fres)
	{
	  if(processID == 0)
	    printf("error reading weight file probably encountered a non-integer weight value\n");
	  errorExit(-1);
	}
      wv[i] = w;
      i++;
    }

  if(i != rdta->sites)
    {
      if(processID == 0)
	printf("number %d of weights not equal to number %d of alignment columns\n", i, rdta->sites);
      errorExit(-1);
    }

  for(i = 1; i <= rdta->sites; i++)
    rdta->wgt[i] = wv[i - 1];

  fclose(weightFile);
  free(wv);
}
static hashNumberType  hashString(char *p, hashNumberType tableSize)
{
  hashNumberType h = 0;
  
  for(; *p; p++)
    h = 31 * h + *p;
  
  return (h % tableSize);
}
static void addword(char *s, stringHashtable *h, int nodeNumber)
{
  hashNumberType position = hashString(s, h->tableSize);
  stringEntry *p = h->table[position];
  
  for(; p!= NULL; p = p->next)
    {
      if(strcmp(s, p->word) == 0)		 
	return;	  	
    }

  p = (stringEntry *)malloc(sizeof(stringEntry));

  assert(p);
  
  p->nodeNumber = nodeNumber;
  p->word = (char *)malloc((strlen(s) + 1) * sizeof(char));

  strcpy(p->word, s);
  
  p->next =  h->table[position];
  
  h->table[position] = p;
}


static stringHashtable *initStringHashTable(hashNumberType n)
{
  /* 
     init with primes 
  */
    
  static const hashNumberType initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
					     196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
					     50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
 

  /* init with powers of two

  static const  hashNumberType initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
					      32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
					      4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
					      268435456, 536870912, 1073741824, 2147483648U};
  */
  
  stringHashtable *h = (stringHashtable*)malloc(sizeof(stringHashtable));
  
  hashNumberType
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]),
    maxSize = (hashNumberType)-1;    

  assert(n <= maxSize);

  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;

  assert(i < primeTableLength);

  tableSize = initTable[i];  

  h->table = (stringEntry**)calloc(tableSize, sizeof(stringEntry*));
  h->tableSize = tableSize;    

  return h;
}



static void getinput(analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
  int i;

  INFILE = myfopen(seq_file, "rb");
  
  getnums(rdta);
  
     
  /*myBinFwrite(&(rdta->sites), sizeof(int), 1);
  myBinFwrite(&(rdta->numsp), sizeof(int), 1);  

  printf("%d %d\n", rdta->sites, rdta->numsp);*/
    

  tr->mxtips            = rdta->numsp;
  
  
  rdta->wgt             = (int *)    malloc((rdta->sites + 1) * sizeof(int));
  cdta->alias           = (int *)    malloc((rdta->sites + 1) * sizeof(int));
  cdta->aliaswgt        = (int *)    malloc((rdta->sites + 1) * sizeof(int)); 
  tr->model             = (int *)    calloc((rdta->sites + 1), sizeof(int));
  tr->initialDataVector  = (int *)    malloc((rdta->sites + 1) * sizeof(int));
  tr->extendedDataVector = (int *)    malloc((rdta->sites + 1) * sizeof(int));         
  
  if(!adef->useWeightFile)
    {
      for (i = 1; i <= rdta->sites; i++)
	rdta->wgt[i] = 1;
    }
  else
    {
      assert(!adef->useSecondaryStructure);
      inputweights(rdta);
    }

  
  if(adef->useMultipleModel)
    {
      int ref;
      
      parsePartitions(adef, rdta, tr);
      
      for(i = 1; i <= rdta->sites; i++)
	{
	  ref = tr->model[i];
	  tr->initialDataVector[i] = tr->initialPartitionData[ref].dataType;
	}
    }
  else
    {
      int dataType = -1;
	              
      tr->initialPartitionData  = (pInfo*)malloc(sizeof(pInfo));
      tr->initialPartitionData[0].partitionName = (char*)malloc(128 * sizeof(char));
      strcpy(tr->initialPartitionData[0].partitionName, "No Name Provided");
      
      tr->initialPartitionData[0].protModels = adef->proteinMatrix;
      tr->initialPartitionData[0].protFreqs  = adef->protEmpiricalFreqs;
      
      
      tr->NumberOfModels = 1;
      
     
      
      if(adef->model == M_PROTCAT || adef->model == M_PROTGAMMA)
	dataType = AA_DATA;
      if(adef->model == M_GTRCAT || adef->model == M_GTRGAMMA)
	dataType = DNA_DATA;
      if(adef->model == M_BINCAT || adef->model == M_BINGAMMA)
	dataType = BINARY_DATA;
      if(adef->model == M_32CAT || adef->model == M_32GAMMA)
	dataType = GENERIC_32;
      if(adef->model == M_64CAT || adef->model == M_64GAMMA)
	dataType = GENERIC_64;
      
      
      
      assert(dataType == BINARY_DATA || dataType == DNA_DATA || dataType == AA_DATA || 
	     dataType == GENERIC_32  || dataType == GENERIC_64);
      
      tr->initialPartitionData[0].dataType = dataType;
      
      for(i = 0; i <= rdta->sites; i++)
	{
	  tr->initialDataVector[i] = dataType;
	  tr->model[i]      = 0;
	}
    }
  
  if(adef->useSecondaryStructure)
    {
      memcpy(tr->extendedDataVector, tr->initialDataVector, (rdta->sites + 1) * sizeof(int));
      
      tr->extendedPartitionData =(pInfo*)malloc(sizeof(pInfo) * tr->NumberOfModels);
      
      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  tr->extendedPartitionData[i].partitionName = (char*)malloc((strlen(tr->initialPartitionData[i].partitionName) + 1) * sizeof(char));
	  strcpy(tr->extendedPartitionData[i].partitionName, tr->initialPartitionData[i].partitionName);
	  tr->extendedPartitionData[i].dataType   = tr->initialPartitionData[i].dataType;
	  
	  tr->extendedPartitionData[i].protModels = tr->initialPartitionData[i].protModels;
	  tr->extendedPartitionData[i].protFreqs  = tr->initialPartitionData[i].protFreqs;
	}
      
      parseSecondaryStructure(tr, adef, rdta->sites);
      
      tr->dataVector    = tr->extendedDataVector;
      tr->partitionData = tr->extendedPartitionData;
    }
  else
    {
      tr->dataVector    = tr->initialDataVector;
      tr->partitionData = tr->initialPartitionData;
    }
  
 
  
  getyspace(rdta);

  setupTree(tr, adef);

      
	
  if(!getdata(adef, rdta, tr))
    {
      printf("Problem reading alignment file \n");
      errorExit(1);
    }
      
  tr->nameHash = initStringHashTable(10 * tr->mxtips);
  for(i = 1; i <= tr->mxtips; i++)
    addword(tr->nameList[i], tr->nameHash, i);
      
  fclose(INFILE);
}



static unsigned char buildStates(int secModel, unsigned char v1, unsigned char v2)
{
  unsigned char new = 0;

  switch(secModel)
    {
    case SECONDARY_DATA:
      new = v1;
      new = new << 4;
      new = new | v2;
      break;
    case SECONDARY_DATA_6:
      {
	int
	  meaningDNA[256],
	  i;

	const unsigned char
	  allowedStates[6][2] = {{'A','T'}, {'C', 'G'}, {'G', 'C'}, {'G','T'}, {'T', 'A'}, {'T', 'G'}};

	const unsigned char
	  finalBinaryStates[6] = {1, 2, 4, 8, 16, 32};

	unsigned char
	  intermediateBinaryStates[6];

	int length = 6;

	for(i = 0; i < 256; i++)
	  meaningDNA[i] = -1;

	meaningDNA['A'] =  1;
	meaningDNA['B'] = 14;
	meaningDNA['C'] =  2;
	meaningDNA['D'] = 13;
	meaningDNA['G'] =  4;
	meaningDNA['H'] = 11;
	meaningDNA['K'] = 12;
	meaningDNA['M'] =  3;
	meaningDNA['N'] = 15;
	meaningDNA['O'] = 15;
	meaningDNA['R'] =  5;
	meaningDNA['S'] =  6;
	meaningDNA['T'] =  8;
	meaningDNA['U'] =  8;
	meaningDNA['V'] =  7;
	meaningDNA['W'] =  9;
	meaningDNA['X'] = 15;
	meaningDNA['Y'] = 10;
	meaningDNA['-'] = 15;
	meaningDNA['?'] = 15;

	for(i = 0; i < length; i++)
	  {
	    unsigned char n1 = meaningDNA[allowedStates[i][0]];
	    unsigned char n2 = meaningDNA[allowedStates[i][1]];

	    new = n1;
	    new = new << 4;
	    new = new | n2;

	    intermediateBinaryStates[i] = new;
	  }

	new = v1;
	new = new << 4;
	new = new | v2;

	for(i = 0; i < length; i++)
	  {
	    if(new == intermediateBinaryStates[i])
	      break;
	  }
	if(i < length)
	  new = finalBinaryStates[i];
	else
	  {
	    new = 0;
	    for(i = 0; i < length; i++)
	      {
		if(v1 & meaningDNA[allowedStates[i][0]])
		  {
		    /*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
		    new |= finalBinaryStates[i];
		  }
		if(v2 & meaningDNA[allowedStates[i][1]])
		  {
		    /*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
		    new |= finalBinaryStates[i];
		  }
	      }
	  }	
      }
      break;
    case SECONDARY_DATA_7:
      {
	int
	  meaningDNA[256],
	  i;

	const unsigned char
	  allowedStates[6][2] = {{'A','T'}, {'C', 'G'}, {'G', 'C'}, {'G','T'}, {'T', 'A'}, {'T', 'G'}};

	const unsigned char
	  finalBinaryStates[7] = {1, 2, 4, 8, 16, 32, 64};

	unsigned char
	  intermediateBinaryStates[7];

	for(i = 0; i < 256; i++)
	  meaningDNA[i] = -1;

	meaningDNA['A'] =  1;
	meaningDNA['B'] = 14;
	meaningDNA['C'] =  2;
	meaningDNA['D'] = 13;
	meaningDNA['G'] =  4;
	meaningDNA['H'] = 11;
	meaningDNA['K'] = 12;
	meaningDNA['M'] =  3;
	meaningDNA['N'] = 15;
	meaningDNA['O'] = 15;
	meaningDNA['R'] =  5;
	meaningDNA['S'] =  6;
	meaningDNA['T'] =  8;
	meaningDNA['U'] =  8;
	meaningDNA['V'] =  7;
	meaningDNA['W'] =  9;
	meaningDNA['X'] = 15;
	meaningDNA['Y'] = 10;
	meaningDNA['-'] = 15;
	meaningDNA['?'] = 15;
	

	for(i = 0; i < 6; i++)
	  {
	    unsigned char n1 = meaningDNA[allowedStates[i][0]];
	    unsigned char n2 = meaningDNA[allowedStates[i][1]];

	    new = n1;
	    new = new << 4;
	    new = new | n2;

	    intermediateBinaryStates[i] = new;
	  }

	new = v1;
	new = new << 4;
	new = new | v2;

	for(i = 0; i < 6; i++)
	  {
	    /* exact match */
	    if(new == intermediateBinaryStates[i])
	      break;
	  }
	if(i < 6)
	  new = finalBinaryStates[i];
	else
	  {
	    /* distinguish between exact mismatches and partial mismatches */

	    for(i = 0; i < 6; i++)
	      if((v1 & meaningDNA[allowedStates[i][0]]) && (v2 & meaningDNA[allowedStates[i][1]]))
		break;
	    if(i < 6)
	      {
		/* printf("partial mismatch\n"); */

		new = 0;
		for(i = 0; i < 6; i++)
		  {
		    if((v1 & meaningDNA[allowedStates[i][0]]) && (v2 & meaningDNA[allowedStates[i][1]]))
		      {
			/*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
			new |= finalBinaryStates[i];
		      }
		    else
		      new |=  finalBinaryStates[6];
		  }
	      }
	    else
	      new = finalBinaryStates[6];
	  }	
      }
      break;
    default:
      assert(0);
    }

  return new;

}



static void adaptRdataToSecondary(tree *tr, rawdata *rdta)
{
  int *alias = (int*)calloc(rdta->sites, sizeof(int));
  int i, j, realPosition;  

  for(i = 0; i < rdta->sites; i++)
    alias[i] = -1;

  for(i = 0, realPosition = 0; i < rdta->sites; i++)
    {
      int partner = tr->secondaryStructurePairs[i];
      if(partner != -1)
	{
	  assert(tr->dataVector[i+1] == SECONDARY_DATA || tr->dataVector[i+1] == SECONDARY_DATA_6 || tr->dataVector[i+1] == SECONDARY_DATA_7);

	  if(i < partner)
	    {
	      for(j = 1; j <= rdta->numsp; j++)
		{
		  unsigned char v1 = rdta->y[j][i+1];
		  unsigned char v2 = rdta->y[j][partner+1];

		  assert(i+1 < partner+1);

		  rdta->y[j][i+1] = buildStates(tr->dataVector[i+1], v1, v2);
		}
	      alias[realPosition] = i;
	      realPosition++;
	    }
	}
      else
	{
	  alias[realPosition] = i;
	  realPosition++;
	}
    }

  assert(rdta->sites - realPosition == tr->numberOfSecondaryColumns / 2);

  rdta->sites = realPosition;

  for(i = 0; i < rdta->sites; i++)
    {
      assert(alias[i] != -1);
      tr->model[i+1]    = tr->model[alias[i]+1];
      tr->dataVector[i+1] = tr->dataVector[alias[i]+1];
      rdta->wgt[i+1] =  rdta->wgt[alias[i]+1];

      for(j = 1; j <= rdta->numsp; j++)
	rdta->y[j][i+1] = rdta->y[j][alias[i]+1];
    }

  free(alias);
}

static void sitesort(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  gap, i, j, jj, jg, k, n, nsp;
  int  
    *index, 
    *category = (int*)NULL;

  boolean  flip, tied;
  unsigned char  **data;

  if(adef->useSecondaryStructure)
    {
      assert(tr->NumberOfModels > 1 && adef->useMultipleModel);

      adaptRdataToSecondary(tr, rdta);
    }

  if(adef->useMultipleModel)    
    category      = tr->model;
  

  index    = cdta->alias;
  data     = rdta->y;
  n        = rdta->sites;
  nsp      = rdta->numsp;
  index[0] = -1;


  if(adef->compressPatterns)
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
		  if(adef->useMultipleModel)
		    {		     		      
		      assert(category[jj] != -1 &&
			     category[jg] != -1);
		     
		      flip = (category[jj] > category[jg]);
		      tied = (category[jj] == category[jg]);		     

		    }
		  else
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


static void sitecombcrunch (rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  i, sitei, j, sitej, k;
  boolean  tied;
  int 
    *aliasModel = (int*)NULL,
    *aliasSuperModel = (int*)NULL;

  if(adef->useMultipleModel)
    {
      aliasSuperModel = (int*)malloc(sizeof(int) * (rdta->sites + 1));
      aliasModel      = (int*)malloc(sizeof(int) * (rdta->sites + 1));
    } 

  i = 0;
  cdta->alias[0]    = cdta->alias[1];
  cdta->aliaswgt[0] = 0;

  if(adef->mode == PER_SITE_LL)
    {
      int i;

      assert(0);

      tr->patternPosition = (int*)malloc(sizeof(int) * rdta->sites);
      tr->columnPosition  = (int*)malloc(sizeof(int) * rdta->sites);

      for(i = 0; i < rdta->sites; i++)
	{
	  tr->patternPosition[i] = -1;
	  tr->columnPosition[i]  = -1;
	}
    }

  

  i = 0;
  for (j = 1; j <= rdta->sites; j++)
    {
      sitei = cdta->alias[i];
      sitej = cdta->alias[j];
      if(!adef->compressPatterns)
	tied = 0;
      else
	{
	  if(adef->useMultipleModel)
	    {	     
	      tied = (tr->model[sitei] == tr->model[sitej]);
	      if(tied)
		assert(tr->dataVector[sitei] == tr->dataVector[sitej]);
	    }
	  else
	    tied = 1;
	}

      for (k = 1; tied && (k <= rdta->numsp); k++)
	tied = (rdta->y[k][sitei] == rdta->y[k][sitej]);

      if (tied)
	{
	  if(adef->mode == PER_SITE_LL)
	    {
	      tr->patternPosition[j - 1] = i;
	      tr->columnPosition[j - 1] = sitej;
	      /*printf("Pattern %d from column %d also at site %d\n", i, sitei, sitej);*/
	    }


	  cdta->aliaswgt[i] += rdta->wgt[sitej];
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
      else
	{
	  if (cdta->aliaswgt[i] > 0) i++;

	  if(adef->mode == PER_SITE_LL)
	    {
	      tr->patternPosition[j - 1] = i;
	      tr->columnPosition[j - 1] = sitej;
	      /*printf("Pattern %d is from cloumn %d\n", i, sitej);*/
	    }

	  cdta->aliaswgt[i] = rdta->wgt[sitej];
	  cdta->alias[i] = sitej;
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
    }

  cdta->endsite = i;
  if (cdta->aliaswgt[i] > 0) cdta->endsite++;

  if(adef->mode == PER_SITE_LL)
    {
      assert(0);

      for(i = 0; i < rdta->sites; i++)
	{
	  int p  = tr->patternPosition[i];
	  int c  = tr->columnPosition[i];

	  assert(p >= 0 && p < cdta->endsite);
	  assert(c >= 1 && c <= rdta->sites);
	}
    }


  if(adef->useMultipleModel)
    {
      for(i = 0; i <= rdta->sites; i++)
	{
	  tr->model[i]      = aliasModel[i];
	  tr->dataVector[i] = aliasSuperModel[i];
	}
    }

  if(adef->useMultipleModel)
    {
      free(aliasModel);
      free(aliasSuperModel);
    }     
}


static boolean makeweights (analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
  int  i;

 
    
  for (i = 1; i <= rdta->sites; i++)
    cdta->alias[i] = i;

  sitesort(rdta, cdta, tr, adef);
  sitecombcrunch(rdta, cdta, tr, adef);
  
      
  /*myBinFwrite(cdta->alias,    sizeof(int), (rdta->sites + 1));
  myBinFwrite(cdta->aliaswgt, sizeof(int), (rdta->sites + 1));
  myBinFwrite(tr->model,      sizeof(int), (rdta->sites + 1));
  myBinFwrite(tr->dataVector, sizeof(int), (rdta->sites + 1));
  myBinFwrite(&(cdta->endsite), sizeof(int), 1);   */

  return TRUE;
}




static boolean makevalues(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  
    i, 
    j, 
    model, 
    modelCounter;

  unsigned char
    *y    = (unsigned char *)malloc(((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char));
  

  /*

  printf("compressed data Assigning %Zu bytes\n", ((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char));

  */
  
  
    {
      for (i = 1; i <= rdta->numsp; i++)
	for (j = 0; j < cdta->endsite; j++)   
	  y[(((size_t)(i - 1)) * ((size_t)cdta->endsite)) + j] = rdta->y[i][cdta->alias[j]];
      
      /*
	printf("Free on raw data\n");
      */

      free(rdta->y0);
      free(rdta->y);
      
   
      /*myBinFwrite(y, sizeof(unsigned char), ((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char));*/
    }

  rdta->y0 = y;
 
  if(!adef->useMultipleModel)
    tr->NumberOfModels = 1;

  if(adef->useMultipleModel)
    {
      tr->partitionData[0].lower = 0;

      model        = tr->model[0];
      modelCounter = 0;
     
      i            = 1;

      while(i <  cdta->endsite)
	{
	  if(tr->model[i] != model)
	    {
	      tr->partitionData[modelCounter].upper     = i;
	      tr->partitionData[modelCounter + 1].lower = i;

	      model = tr->model[i];	     
	      modelCounter++;
	    }
	  i++;
	}


      tr->partitionData[tr->NumberOfModels - 1].upper = cdta->endsite;      
    
      for(i = 0; i < tr->NumberOfModels; i++)		  
	tr->partitionData[i].width      = tr->partitionData[i].upper -  tr->partitionData[i].lower;
	 
      model        = tr->model[0];
      modelCounter = 0;
      tr->model[0] = modelCounter;
      i            = 1;
	
      while(i < cdta->endsite)
	{	 
	  if(tr->model[i] != model)
	    {
	      model = tr->model[i];
	      modelCounter++;
	      tr->model[i] = modelCounter;
	    }
	  else
	    tr->model[i] = modelCounter;
	  i++;
	}      
    }
  else
    {
      tr->partitionData[0].lower = 0;
      tr->partitionData[0].upper = cdta->endsite;
      tr->partitionData[0].width =  tr->partitionData[0].upper -  tr->partitionData[0].lower;
    }

  tr->rdta       = rdta;
  tr->cdta       = cdta; 

  tr->originalCrunchedLength = tr->cdta->endsite;
    
  for(i = 0; i < rdta->numsp; i++)
    tr->yVector[i + 1] = &(rdta->y0[((size_t)tr->originalCrunchedLength) * ((size_t)i)]);

  return TRUE;
}



static void initAdef(analdef *adef)
{  
  adef->useSecondaryStructure  = FALSE;
  adef->bootstrapBranchLengths = FALSE;
  adef->model                  = M_GTRCAT;
  adef->max_rearrange          = 21;
  adef->stepwidth              = 5;
  adef->initial                = adef->bestTrav = 10;
  adef->initialSet             = FALSE;
  adef->restart                = FALSE;
  adef->mode                   = BIG_RAPID_MODE;
  adef->categories             = 25;
  adef->boot                   = 0;
  adef->rapidBoot              = 0;
  adef->useWeightFile          = FALSE;
  adef->checkpoints            = 0;
  adef->startingTreeOnly       = 0;
  adef->multipleRuns           = 1;
  adef->useMultipleModel       = FALSE;
  adef->likelihoodEpsilon      = 0.1;
  adef->constraint             = FALSE;
  adef->grouping               = FALSE;
  adef->randomStartingTree     = FALSE;
  adef->parsimonySeed          = 0;
  adef->proteinMatrix          = JTT;
  adef->protEmpiricalFreqs     = 0;  
  adef->useInvariant           = FALSE;
  adef->permuteTreeoptimize    = FALSE;
  adef->useInvariant           = FALSE;
  adef->allInOne               = FALSE;
  adef->likelihoodTest         = FALSE;
  adef->perGeneBranchLengths   = FALSE;
  adef->generateBS             = FALSE;
  adef->bootStopping           = FALSE;
  adef->gapyness               = 0.0;
  adef->similarityFilterMode   = 0;
  adef->useExcludeFile         = FALSE;
  adef->userProteinModel       = FALSE;
  adef->externalAAMatrix       = (double*)NULL;
  adef->computeELW             = FALSE;
  adef->computeDistance        = FALSE;
  adef->thoroughInsertion      = FALSE;
  adef->compressPatterns       = TRUE; 
  adef->readTaxaOnly           = FALSE;
  adef->meshSearch             = 0;
  adef->useCheckpoint          = FALSE;
  adef->leaveDropMode          = FALSE;
  adef->slidingWindowSize      = 100;
#ifdef _BAYESIAN 
  adef->bayesian               = FALSE;
#endif

}




static int dataExists(char *model, analdef *adef)
{
  int i;
  char thisModel[1024];

  /********** BINARY ********************/

   if(strcmp(model, "BIN\0") == 0)
    {
      adef->model = M_BINGAMMA;      
      return 1;
    }

  

  /*********** 32 state ****************************/

  if(strcmp(model, "MULTI\0") == 0)
    {
      adef->model = M_32GAMMA;     
      return 1;
    }
  

  /*********** 64 state ****************************/

  if(strcmp(model, "CODON\0") == 0)
    {
      adef->model = M_64GAMMA;     
      return 1;
    }

  

  /*********** DNA **********************/

  if(strcmp(model, "DNA\0") == 0)
    {
      adef->model = M_GTRGAMMA;     
      return 1;
    }

 



  /*************** AA GTR ********************/

  /* TODO empirical FREQS */

  if(strcmp(model, "PROT\0") == 0)
    {
      adef->model = M_PROTGAMMA;     
      return 1;
    } 

 



  return 0;
}

/*********************************************************************************************/

static void printVersionInfo(void)
{
  printf("\n\nThis is %s version %s released by Alexandros Stamatakis on %s.\n\n",  programName, programVersion, programDate); 
}

static void printREADME(void)
{
  printVersionInfo();
  printf("\n");  
  printf("\nTo report bugs send an email to raxml@h-its.org\n");
  printf("Please send me all input files, the exact invocation, details of the HW and operating system,\n");
  printf("as well as all error messages printed to screen.\n\n\n");

  printf("parser\n");
  printf("      -s sequenceFileName\n");
  printf("      -n outputFileName\n");
  printf("      -m substitutionModel\n");
  printf("      [-c]\n");
  printf("      [-q]\n");
  printf("      [-h]\n");
  printf("\n"); 
  printf("      -m Model of  Nucleotide or Amino Acid Substitution:\n");
  printf("\n"); 
  printf("              For DNA data use:    DNA\n");	
  printf("              For AA data use:     PROT\n");			   
  printf("\n"); 
  printf("      -c      disable site pattern compression\n");
  printf("\n");
  printf("      -q      Specify the file name which contains the assignment of models to alignment\n");
  printf("              partitions for multiple models of substitution. For the syntax of this file\n");
  printf("              please consult the manual.\n");  
  printf("\n");
  printf("      -h      Display this help message.\n");
  printf("\n");
  printf("\n\n\n\n");

}

static int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg)
{
  static int sp = 1;
  register int c;
  register char *cp;

  if(sp == 1)
    {
      if(*optind >= argc || argv[*optind][0] != '-' || argv[*optind][1] == '\0')
        {
	  return -1;
        }
    }
  else
    {
      if(strcmp(argv[*optind], "--") == 0)
	{
	  *optind =  *optind + 1;
	  return -1;
	}
    }

  c = argv[*optind][sp];
  if(c == ':' || (cp=strchr(opts, c)) == 0)
    {
      if(argv[*optind][++sp] == '\0')
	{
	  *optind =  *optind + 1;
	  sp = 1;
	}
      printf("\n Error: illegal option -- %c\n\n", c);
      return('?');
    }
  if(*++cp == ':')
    {
      if(argv[*optind][sp+1] != '\0')
	{
	  *optarg = &argv[*optind][sp+1];
	  *optind =  *optind + 1;
	}
      else
	{
	  *optind =  *optind + 1;
	  if(*optind >= argc)
	    {
	      if ( c != 'h')	
                {
	          sp = 1;
                  printf("\n Error: option -- %c requires an argument\n\n", c);
	          return('?');
                }
               else
                  return ( c );
	    }
	  else
	    {
	      *optarg = argv[*optind];
	      *optind =  *optind + 1;
	    }
	}
      sp = 1;
    }
  else
    {
      if(argv[*optind][++sp] == '\0')
	{
	  sp = 1;
	  *optind =  *optind + 1;
	}
      *optarg = 0;
    }
  return(c);
  }


/*********************************************************************************************/









static void analyzeRunId(char id[128])
{
  int i = 0;

  while(id[i] != '\0')
    {    
      if(i >= 128)
	{
	  printf("\n Error: run id after \"-n\" is too long, it has %d characters please use a shorter one\n\n", i);
	  assert(0);
	}
      
      if(id[i] == '/')
	{
	  printf("\n Error character %c not allowed in run ID\n\n", id[i]);
	  assert(0);
	}


      i++;
    }

  if(i == 0)
    {
      printf("\n Error: please provide a string for the run id after \"-n\" \n\n");
      assert(0);
    }

}


static void get_args(int argc, char *argv[], analdef *adef, tree *tr)
{
  boolean
    bad_opt    =FALSE,
    resultDirSet = FALSE;

  char
    resultDir[1024] = "",          
    *optarg,
    model[2048] = "",       
    modelChar;

  double 
    likelihoodEpsilon;
  
  int  
    optind = 1,        
    c,
    nameSet = 0,
    alignmentSet = 0,    
    treeSet = 0,   
    modelSet = 0;


  run_id[0] = 0; 
  seq_file[0] = 0;
  model[0] = 0;
  weightFileName[0] = 0;
  modelFileName[0] = 0;

  /*********** tr inits **************/

#ifdef _USE_PTHREADS
  NumberOfThreads = 0;
#endif
  
 
  tr->bootStopCriterion = -1;
  tr->wcThreshold = 0.03;
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->catOnly = FALSE;
 
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->useGappedImplementation = FALSE;
  tr->saveMemory = FALSE;
  


  
  /********* tr inits end*************/


    while( !bad_opt && ( ( c = mygetopt(argc,argv,"q:s:n:m:hc", &optind, &optarg ) ) != -1 ) )
    {
    switch(c)
      {                
      case 'c':
	adef->compressPatterns = FALSE;
	break;
      case 'h':
        printREADME();
	errorExit(-1);
        break;                 
      case 'q':
	strcpy(modelFileName,optarg);
	adef->useMultipleModel = TRUE;
        break;                 
      case 'n':
        strcpy(run_id,optarg);
	analyzeRunId(run_id);
	nameSet = 1;
        break;     
      case 's':
	strcpy(seq_file, optarg);
	alignmentSet = 1;
	break;
      case 'm':
	strcpy(model,optarg);
	if(dataExists(model, adef) == 0)
	  {
	    printf("\n Error: model %s does not exist\n\n", model);               
	    errorExit(-1);
	  }
	else
	  modelSet = 1;
	break; 
      default:
	errorExit(-1);
    }
  }  

  if(!adef->useMultipleModel && !modelSet)
    {
      if(processID == 0)
        {
          printREADME();	    
	  printf("\n Error, you must specify a data type for unpartitioned alignment with the \"-m\" option\n\n");
        }
      errorExit(-1);
    }

  if(!nameSet)
    {
      if(processID == 0)
        {
          printREADME();	    
	  printf("\n Error: please specify a name for this run with -n\n\n");
        }
      errorExit(-1);
    }


  if(!alignmentSet)
    {
      if(processID == 0)
        {
          printREADME();	    
	  printf("\n Error: please specify an alignment for this run with -s\n\n");
        } 
      errorExit(-1);
    }
  
  
   strcat(infoFileName,         "RAxML_info."); 
   strcat(infoFileName,         run_id);
  
   if(processID == 0)
     {
       int infoFileExists = 0;
       
       infoFileExists = filexists(infoFileName);
       
       if(infoFileExists)
	 {
	   printf("\n Error: output files with the run ID <%s> already exist... exiting\n\n", run_id);
	   exit(-1);
	 }
     }

  strcat(byteFileName, run_id);
  strcat(byteFileName, ".binary");
  
  if(filexists(byteFileName))
    {
      printf("\n Error: binary compressed file %s you want to generate already exists... exiting\n\n", byteFileName);
      exit(0);
    }

  byteFile = fopen(byteFileName, "wb");
  if ( !byteFile )  
    printf("%s\n", byteFileName);

  return;
}




void errorExit(int e)
{

#ifdef _WAYNE_MPI
  MPI_Finalize();
#endif

  exit(e);

}





 




/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/





void getDataTypeString(tree *tr, int model, char typeOfData[1024])
{
  switch(tr->partitionData[model].dataType)
    {
    case AA_DATA:
      strcpy(typeOfData,"AA");
      break;
    case DNA_DATA:
      strcpy(typeOfData,"DNA");
      break;
    case BINARY_DATA:
      strcpy(typeOfData,"BINARY/MORPHOLOGICAL");
      break;
    case SECONDARY_DATA:
      strcpy(typeOfData,"SECONDARY 16 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_6:
      strcpy(typeOfData,"SECONDARY 6 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_7:
      strcpy(typeOfData,"SECONDARY 7 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case GENERIC_32:
      strcpy(typeOfData,"Multi-State");
      break;
    case GENERIC_64:
      strcpy(typeOfData,"Codon"); 
      break;
    default:
      assert(0);
    }
}





/************************************************************************************/





  



static int iterated_bitcount(unsigned int n)
{
    int 
      count=0;    
    
    while(n)
      {
        count += n & 0x1u ;    
        n >>= 1 ;
      }
    
    return count;
}

static char bits_in_16bits [0x1u << 16];

static void compute_bits_in_16bits(void)
{
    unsigned int i;    
    
    for (i = 0; i < (0x1u<<16); i++)
        bits_in_16bits[i] = iterated_bitcount(i);
    
    return ;
}

unsigned int precomputed16_bitcount (unsigned int n)
{
  /* works only for 32-bit int*/
    
    return bits_in_16bits [n         & 0xffffu]
        +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}







static void smoothFreqs(const int n, double *pfreqs, double *dst, pInfo *partitionData)
{
  int 
    countScale = 0, 
    l,
    loopCounter = 0;  
  

  /*
    for(l = 0; l < n; l++)
    if(pfreqs[l] < FREQ_MIN)
      countScale++;
  */

  for(l = 0; l < n; l++)
    if(pfreqs[l] == 0.0)
      countScale++;

  if(countScale > 0)
    {	     
      while(countScale > 0)
	{
	  double correction = 0.0;
	  double factor = 1.0;
	  
	  for(l = 0; l < n; l++)
	    {
	      if(pfreqs[l] == 0.0)		  
		correction += FREQ_MIN;		   		  
	      else
		if(pfreqs[l] < FREQ_MIN)		    
		  {
		    correction += (FREQ_MIN - pfreqs[l]);
		    factor -= (FREQ_MIN - pfreqs[l]);
		  }
	    }		      	    	    
	  
	  countScale = 0;
	  
	  for(l = 0; l < n; l++)
	    {		    
	      if(pfreqs[l] >= FREQ_MIN)		      
		pfreqs[l] = pfreqs[l] - (pfreqs[l] * correction * factor);	
	      else
		pfreqs[l] = FREQ_MIN;
	      
	      if(pfreqs[l] < FREQ_MIN)
		countScale++;
	    }
	  assert(loopCounter < 100);
	  loopCounter++;
	}		    
    }

  for(l = 0; l < n; l++)
    dst[l] = pfreqs[l];

  
  if(partitionData->nonGTR)
    {
      int k;

      assert(partitionData->dataType == SECONDARY_DATA_7 || partitionData->dataType == SECONDARY_DATA_6 || partitionData->dataType == SECONDARY_DATA);
       
      for(l = 0; l < n; l++)
	{
	  int count = 1;	
	  
	  for(k = 0; k < n; k++)
	    {
	      if(k != l && partitionData->frequencyGrouping[l] == partitionData->frequencyGrouping[k])
		{
		  count++;
		  dst[l] += pfreqs[k];
		}
	    }
	  dst[l] /= ((double)count);
	}            
     }  
}
	    

static void genericBaseFrequencies(tree *tr, const int numFreqs, rawdata *rdta, cruncheddata *cdta, int lower, int upper, int model, boolean smoothFrequencies,
				   const unsigned int *bitMask)
{
  double 
    wj, 
    acc,
    pfreqs[64], 
    sumf[64],   
    temp[64];
 
  int     
    i, 
    j, 
    k, 
    l;

  unsigned char  *yptr;  
	  
  for(l = 0; l < numFreqs; l++)	    
    pfreqs[l] = 1.0 / ((double)numFreqs);
	  
  for (k = 1; k <= 8; k++) 
    {	     	   	    	      			
      for(l = 0; l < numFreqs; l++)
	sumf[l] = 0.0;
	      
      for (i = 0; i < rdta->numsp; i++) 
	{		 
	  yptr =  &(rdta->y0[((size_t)i) * ((size_t)tr->originalCrunchedLength)]);
	  
	  for(j = lower; j < upper; j++) 
	    {
	      unsigned int code = bitMask[yptr[j]];
	      assert(code >= 1);
	      
	      for(l = 0; l < numFreqs; l++)
		{
		  if((code >> l) & 1)
		    temp[l] = pfreqs[l];
		  else
		    temp[l] = 0.0;
		}		      	      
	      
	      for(l = 0, acc = 0.0; l < numFreqs; l++)
		{
		  if(temp[l] != 0.0)
		    acc += temp[l];
		}
	      
	      wj = ((double)cdta->aliaswgt[j]) / acc;
	      
	      for(l = 0; l < numFreqs; l++)
		{
		  if(temp[l] != 0.0)		    
		    sumf[l] += wj * temp[l];			     				   			     		   
		}
	    }
	}	    	      
      
      for(l = 0, acc = 0.0; l < numFreqs; l++)
	{
	  if(sumf[l] != 0.0)
	    acc += sumf[l];
	}
	      
      for(l = 0; l < numFreqs; l++)
	pfreqs[l] = sumf[l] / acc;	     
    }
  
  if(smoothFrequencies)         
    smoothFreqs(numFreqs, pfreqs,  tr->partitionData[model].frequencies, &(tr->partitionData[model]));	   
  else    
    {
      boolean 
	zeroFreq = FALSE;

      char 
	typeOfData[1024];

      getDataTypeString(tr, model, typeOfData);  

      for(l = 0; l < numFreqs; l++)
	{
	  if(pfreqs[l] == 0.0)
	    {
	      printBothOpen("Empirical base frequency for state number %d is equal to zero in %s data partition %s\n", l, typeOfData, tr->partitionData[model].partitionName);
	      printBothOpen("Since this is probably not what you want to do, RAxML will soon exit.\n\n");
	      zeroFreq = TRUE;
	    }
	}

      if(zeroFreq)
	exit(-1);

      for(l = 0; l < numFreqs; l++)
	{
	  assert(pfreqs[l] > 0.0);
	  tr->partitionData[model].frequencies[l] = pfreqs[l];
	}     
    }  
 
}







static void baseFrequenciesGTR(rawdata *rdta, cruncheddata *cdta, tree *tr)
{  
  int 
    model,
    lower,
    upper,
    states;

  for(model = 0; model < tr->NumberOfModels; model++)
    {      
      lower = tr->partitionData[model].lower;
      upper = tr->partitionData[model].upper;	  	 
      states = tr->partitionData[model].states;
	
      switch(tr->partitionData[model].dataType)
	{
	case GENERIC_32:
	  switch(tr->multiStateModel)
	    {
	    case ORDERED_MULTI_STATE:
	    case MK_MULTI_STATE:	   
	      {	       
		int i;
		double 
		  freq = 1.0 / (double)states,
		  acc = 0.0;

		for(i = 0; i < states; i++)
		  {
		    acc += freq;
		    tr->partitionData[model].frequencies[i] = freq;
		    /*printf("%f \n", freq);*/
		  }
		/*printf("Frequency Deviation: %1.60f\n", acc);*/
	      }
	      break;
	     case GTR_MULTI_STATE:
	      genericBaseFrequencies(tr, states, rdta, cdta, lower, upper, model, TRUE,
				     bitVector32);
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case GENERIC_64:	 
	  assert(0);
	  break;
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case SECONDARY_DATA:
	case AA_DATA:
	case DNA_DATA:
	case BINARY_DATA:	  
	  genericBaseFrequencies(tr, states, rdta, cdta, lower, upper, model, 
				 getSmoothFreqs(tr->partitionData[model].dataType),
				 getBitVector(tr->partitionData[model].dataType));	  	 
	  break;	
	default:
	  assert(0);     
	}      
    }
  
  return;
}



int main (int argc, char *argv[])
{
  int model;

  rawdata      *rdta;
  cruncheddata *cdta;
  tree         *tr;
  analdef      *adef;
  
  /* get the start time */

  masterTime = gettime();

  /* get some memory for the basic data structures */

  adef = (analdef *)malloc(sizeof(analdef));
  rdta = (rawdata *)malloc(sizeof(rawdata));
  cdta = (cruncheddata *)malloc(sizeof(cruncheddata));
  tr   = (tree *)malloc(sizeof(tree));


  /* the initialization below is required for the hash tables that are used */

  compute_bits_in_16bits();

  /* initialize the analysis parameters in struct adef to default values */

  initAdef(adef);

  /* parse command line arguments: this has a side effect on tr struct and adef struct variables */

  get_args(argc,argv, adef, tr); 
            
  /* parse the phylip file: this should probably be re-done, perhaps using the relatively flexible parser 
     written in C++ by Marc Holder */
  
  getinput(adef, rdta, cdta, tr);  

  printBothOpen("Pattern compression: %s\n", (adef->compressPatterns)?"ON":"OFF");

  makeweights(adef, rdta, cdta, tr);         
      
  makevalues(rdta, cdta, tr, adef);                 
                   
                  
  for(model = 0; model < tr->NumberOfModels; model++)
    {	
      int 
	states = -1,
	maxTipStates = getUndetermined(tr->partitionData[model].dataType) + 1;  	      
      
      const 
	partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model]));
      
      tr->partitionData[model].frequencies       = (double*)malloc(pl->frequenciesLength * sizeof(double)); 
            
      switch(tr->partitionData[model].dataType)
	{
	case DNA_DATA:
	case AA_DATA:	
	  states = getStates(tr->partitionData[model].dataType);	 
	  break;	
	default:
	  assert(0);
	}

      tr->partitionData[model].states       = states;
      tr->partitionData[model].maxTipStates = maxTipStates;
    }   
  
  baseFrequenciesGTR(tr->rdta, tr->cdta, tr); 
  
 

  /*for(model = 0; model < tr->NumberOfModels; model++)	    	    
    {
      int i;

     

      myBinFwrite(tr->partitionData[model].frequencies, sizeof(double), tr->partitionData[model].states);	      	   
      for(i = 0; i < tr->partitionData[model].states; i++)
	printf("%f ", tr->partitionData[model].frequencies[i]);
      printf("\n");
      }*/
  

  {
    size_t 
      i,
      model;
    
    unsigned char *y;
    
    myBinFwrite(&(tr->mxtips),                 sizeof(int), 1);
    myBinFwrite(&(tr->originalCrunchedLength), sizeof(int), 1);
    myBinFwrite(&(tr->NumberOfModels),         sizeof(int), 1);
    myBinFwrite(&(adef->gapyness),             sizeof(double), 1);
    myBinFwrite(tr->cdta->aliaswgt,               sizeof(int), tr->originalCrunchedLength);	  	  	       	
	
    for(i = 1; i <= tr->mxtips; i++)
      {
	int len = strlen(tr->nameList[i]) + 1;
	myBinFwrite(&len, sizeof(int), 1);
	myBinFwrite(tr->nameList[i], sizeof(char), len);	
      }  
	  	
    for(model = 0; model < (size_t)tr->NumberOfModels; model++)
      {
	int 
	  len;
	
	pInfo 
	  *p = &(tr->partitionData[model]);
	
	myBinFwrite(&(p->states),             sizeof(int), 1);
	myBinFwrite(&(p->maxTipStates),       sizeof(int), 1);
	myBinFwrite(&(p->lower),              sizeof(int), 1);
	myBinFwrite(&(p->upper),              sizeof(int), 1);
	myBinFwrite(&(p->width),              sizeof(int), 1);
	myBinFwrite(&(p->dataType),           sizeof(int), 1);
	myBinFwrite(&(p->protModels),         sizeof(int), 1);
	myBinFwrite(&(p->autoProtModels),     sizeof(int), 1);
	myBinFwrite(&(p->protFreqs),          sizeof(int), 1);
	myBinFwrite(&(p->nonGTR),             sizeof(boolean), 1);
	myBinFwrite(&(p->numberOfCategories), sizeof(int), 1);	 
	
	/* later on if adding secondary structure data
	   
	   int    *symmetryVector;
	   int    *frequencyGrouping;
	*/
	
	len = strlen(p->partitionName) + 1;
	myBinFwrite(&len, sizeof(int), 1);
	myBinFwrite(p->partitionName, sizeof(char), len);	    
	myBinFwrite(tr->partitionData[model].frequencies, sizeof(double), tr->partitionData[model].states);
      }	            
      
    
    myBinFwrite(rdta->y0, sizeof(unsigned char), ((size_t)tr->originalCrunchedLength) * ((size_t)tr->mxtips));          
  }


  fclose(byteFile);  
  
  printBothOpen("\n\nBinary and compressed alignment file written to file %s\n\n", byteFileName);
  printBothOpen("Parsing completed, exiting now ... \n\n");

  return 0;
}


