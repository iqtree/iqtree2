#include <stdio.h>
#include "stack.h"
#include "mem_alloc.h"

int pllStackSize (pllStack ** stack)
{
  pllStack * top;
  int size = 0;
  top = *stack;
 
  while (top)
  {
    ++ size;
    top = top->next;
  }
  
  return (size);
}

int 
pllStackPush (pllStack ** head, void * item)
{
  pllStack * new;
 
  new = (pllStack *) rax_malloc (sizeof (pllStack));
  if (!new) return (0);
 
  new->item = item;
  new->next = *head;
  *head     = new;
 
  return (1);
}

void * pllStackPop (pllStack ** head)
{
  struct item_t * item;
  pllStack * tmp;
  if (!*head) return (NULL);
 
  tmp     = (*head);
  item    = (*head)->item;
  (*head) = (*head)->next;
  rax_free (tmp);
 
  return (item);
}
 
void 
pllStackClear (pllStack ** stack)
{
  while (*stack) pllStackPop (stack);
}

