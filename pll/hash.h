#ifndef __pll_HASH__
#define __pll_HASH__

struct pllHashItem
{
  void * data;
  char * str;
  struct pllHashItem * next;
};

struct pllHashTable
{
  unsigned int size;
  struct pllHashItem ** Items;
};

unsigned int pllHashString (const char * s, unsigned int size);
int pllHashAdd  (struct pllHashTable * hTable, const char * s, void * item);
struct pllHashTable * pllHashInit (unsigned int n);
int pllHashSearch (struct pllHashTable * hTable, char * s, void ** item);
void pllHashDestroy (struct pllHashTable ** hTable, int);
#endif
