#include <stdio.h>
#include <stdlib.h>
#include "mem_alloc.h"

/*  string sorting implementation from:
 *  Bentley J. L., Sedgewick R.: Fast Algorithms for Sorting and Searching 
 *  Strings. In Proceedings of ACM-SIAM Symposium on Discrete Algorithms 
 *  (SODA) 1997.
 */

static void 
vecswap (int i, int j, int n, char ** x, int * oi)
{
  while (n-- > 0)
   {
     PLL_SWAP (x[i], x[j]);
     PLL_SWAP (oi[i], oi[j]);

     ++ i; ++ j;
   }
}

static void 
ssort1 (char ** x, int n, int depth, int * oi)
{
  int           a, b, c, d, r, v;

  if (n <= 1) return;

  a = rand() % n;

  PLL_SWAP (x[0], x[a]);
  PLL_SWAP (oi[0], oi[a]);

  v = x[0][depth];

  a = b = 1;
  c = d = n - 1;

  for (;;)
   {
     while (b <= c && (r = x[b][depth] - v) <= 0)
      {
        if (r == 0)
         {
           PLL_SWAP (x[a], x[b]);
           PLL_SWAP (oi[a], oi[b]);
           ++ a;
         }
        ++ b;
      }
     while (b <= c && (r = x[c][depth] - v) >= 0)
      {
        if (r == 0)
         {
           PLL_SWAP (x[c], x[d]);
           PLL_SWAP (oi[c], oi[d]);
           -- d;
         }
        -- c;
      }
     if (b > c) break;
     PLL_SWAP (x[b], x[c]);
     PLL_SWAP (oi[b], oi[c]);
     ++ b; -- c;
   }
  r = PLL_MIN (a,     b - a);      vecswap (0, b - r, r, x, oi);
  r = PLL_MIN (d - c, n - d - 1);  vecswap (b, n - r, r, x, oi);
  r = b - a; ssort1 (x, r, depth, oi);
  if (x[r][depth] != 0)
   {
     ssort1 (x + r, a + n - d - 1, depth + 1, oi + r);
   }
  r = d - c; ssort1 (x + n - r, r, depth, oi + n - r);
}

int * 
pllssort1main (char ** x, int n)
{
  int * oi;
  int i;

  oi = (int *) rax_malloc (n * sizeof (int));
  for (i = 0; i < n; ++ i)
   {
     oi[i] = i;
   }
  ssort1 (x, n, 0, oi);
  
  return (oi);
}

