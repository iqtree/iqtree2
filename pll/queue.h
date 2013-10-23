#ifndef __pll_QUEUE__
#define __pll_QUEUE__

struct pllQueueItem
{  
  void * item;
  struct pllQueueItem * next;
}; 
   
typedef struct
{  
  struct pllQueueItem * head;
  struct pllQueueItem * tail;
} pllQueue; 

int pllQueueInit (pllQueue ** q);
int pllQueueSize (pllQueue * q);
int pllQueueRemove (pllQueue * q, void ** item);
int pllQueueAppend (pllQueue * q, void * item);
#endif
