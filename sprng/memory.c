#include <stdio.h>
#include <stdlib.h>

#ifdef __STDC__
void *_mymalloc(long size, int line, char *message)
#else
void *_mymalloc(size, line, message)
long size;
int line;
char *message;
#endif
{
  char *temp;

  if(size == 0)
    return NULL;

  temp = (char *) malloc(size);
  
  if(temp == NULL)
    {
      fprintf(stderr,"\nmemory allocation failure in file: %s at line number: %d\n", message, line);
      return NULL;
    }

  return (void *) temp;
}
