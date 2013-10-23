#ifndef __pll_STACK__
#define __pll_STACK__

struct pllStack
{
  void * item;
  struct pllStack * next;
};

typedef struct pllStack pllStack;

void  pllStackClear (pllStack ** stack);
void * pllStackPop (pllStack ** head);
int pllStackPush (pllStack ** head, void * item);
int pllStackSize (pllStack ** stack);

#endif
