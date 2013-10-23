#ifndef __pll_PART__
#define __pll_PART__
#include "queue.h"

typedef struct
{
  int start;
  int end;
  int stride;
} pllPartitionRegion;

typedef struct 
{
  char * partitionName;
  char * partitionModel;
  int protModels;
  int protFreqs;
  int dataType;
  int optimizeBaseFrequencies;
  pllQueue * regionList;
} pllPartitionInfo;
#endif
