#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "memory.h"

#define CHECK_MASK 0x1c		/* Assumes exactly last two bits are 0 ... */
                                /* ... If not, change checklist dimension  */
typedef struct checkidstruct
{
  int *ID;
  struct checkidstruct *next;
} Checklisttype;

Checklisttype checklist[(CHECK_MASK>>2)+1] = {{NULL,NULL},{NULL,NULL},{NULL,NULL},{NULL,NULL}};


#ifdef __STDC__
int *checkID( int *ptr)
#else
int *checkID(ptr)
int *ptr;
#endif
{
  int bucket;
  Checklisttype *next;
  
  if(ptr == NULL)
    return NULL;
  
  bucket = (((uintptr_t) ptr)&CHECK_MASK)>>2;
  next = checklist[bucket].next;
  
  while(next != NULL)
  {
    if(next->ID == ptr)
    {
      return (int *) ptr;
    }
    else
    {
      next = next->next;
    }
    
  }
  
  fprintf(stderr,"ERROR: Invalid generator ID %p\n", (void *) ptr);
  return NULL;
}



#ifdef __STDC__
int *deleteID( int *ptr)
#else
int *deleteID(ptr)
int *ptr;
#endif
{
  int bucket;
  Checklisttype *next, *temp;
  
  
  if(ptr == NULL)
    return NULL;
  
  bucket = (((uintptr_t) ptr)&CHECK_MASK)>>2;
  next = &checklist[bucket];
  
  while(next->next != NULL)
    if(next->next->ID == ptr)
    {
      temp = next->next;
      next->next = next->next->next;
      
      free(temp);
      return (int *) ptr;
    }
    else
    {
      next = next->next;
    }
  

  fprintf(stderr,"ERROR: Invalid generator ID %p\n", (void *) ptr);
  return NULL;
}


#ifdef __STDC__
int *addID( int *ptr)
#else
int *addID(ptr)
int *ptr;
#endif
{
  int bucket;
  /* Checklisttype *next, *temp; (HAS) */
  Checklisttype *temp;
  
  if(ptr == NULL)
    return NULL;
 
  
  bucket = (((uintptr_t) ptr)&CHECK_MASK)>>2;
  
  temp = (Checklisttype *) mymalloc(sizeof(Checklisttype));
  if(temp == NULL)
    return NULL;
  
  temp->ID = (int *) ptr;
  temp->next = checklist[bucket].next;
  checklist[bucket].next = temp;
  
  
  return (int *) ptr;
}



  
  

