#include <stdio.h>
#include "queue.h"
#include "mem_alloc.h"

int
pllQueueInit (pllQueue ** q)
{  
  *q = (pllQueue *) rax_malloc (sizeof (pllQueue));
  if (!*q) return (0);
   
  (*q)->head = NULL;
  (*q)->tail = NULL;
   
  return (1);
}  

int 
pllQueueSize (pllQueue * q)
{  
  int n = 0;
  struct pllQueueItem * elm;
   
  if (!q) return (0);
   
  for (elm = q->head; elm; elm = elm->next) ++n;
   
  return (n);
}  

int
pllQueueRemove (pllQueue * q, void ** item)
{  
  struct pllQueueItem * elm;
   
  if (!q || !q->head) return (0);
   
  elm = q->head;
   
  *item = elm->item;
   
  q->head = q->head->next;
  if (!q->head)  q->tail = NULL;
  rax_free (elm);
   
  return (1);
}  

int 
pllQueueAppend (pllQueue * q, void * item)
{ 
  struct pllQueueItem * qitem;
  if (!q) return (0);
  
  qitem = (struct pllQueueItem *) rax_malloc (sizeof (struct pllQueueItem));
  if (!qitem) return (0);
  
  qitem->item = item;
  qitem->next = NULL;
  
  if (!q->head) 
    q->head = qitem;
  else
    q->tail->next = qitem;
  
  q->tail = qitem;

  return (1);
} 
