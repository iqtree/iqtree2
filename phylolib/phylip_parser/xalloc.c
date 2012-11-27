#include <stdio.h>
#include <stdarg.h>
#include "xalloc.h"


void
xalloc_die (void)
{
  error (EXIT_FAILURE, 0, "%s", "memory exhausted");

  abort();
}

void *
xmalloc (size_t n)
{
  void * p = malloc (n);

  if (!p && n)
     xalloc_die ();

  return (p);
}

void
FREE (int cnt, ...)
{
  va_list ap;
  int i;

  va_start (ap, cnt);
  for (i = 0; i < cnt; ++i)
   {
     free (va_arg (ap, void *));
   }
}
