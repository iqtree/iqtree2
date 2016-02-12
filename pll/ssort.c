/** 
 * PLL (version 1.0.0) a software library for phylogenetic inference
 * Copyright (C) 2013 Tomas Flouri and Alexandros Stamatakis
 *
 * Derived from 
 * RAxML-HPC, a program for sequential and parallel estimation of phylogenetic
 * trees by Alexandros Stamatakis
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For any other enquiries send an Email to Tomas Flouri
 * Tomas.Flouri@h-its.org
 *
 * When publishing work that uses PLL please cite PLL
 * 
 * @file ssort.c
 * Detailed description to appear soon.
 */
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
     PLL_SWAP_PTR (x[i], x[j]);
     PLL_SWAP_INT (oi[i], oi[j]);

     ++ i; ++ j;
   }
}

static void 
ssort1 (char ** x, int n, int depth, int * oi)
{
  int           a, b, c, d, r, v;

  if (n <= 1) return;

  a = rand() % n;

  PLL_SWAP_PTR (x[0], x[a]);
  PLL_SWAP_INT (oi[0], oi[a]);

  v = x[0][depth];

  a = b = 1;
  c = d = n - 1;

  for (;;)
   {
     while (b <= c && (r = x[b][depth] - v) <= 0)
      {
        if (r == 0)
         {
           PLL_SWAP_PTR (x[a], x[b]);
           PLL_SWAP_INT (oi[a], oi[b]);
           ++ a;
         }
        ++ b;
      }
     while (b <= c && (r = x[c][depth] - v) >= 0)
      {
        if (r == 0)
         {
           PLL_SWAP_PTR (x[c], x[d]);
           PLL_SWAP_INT (oi[c], oi[d]);
           -- d;
         }
        -- c;
      }
     if (b > c) break;
     PLL_SWAP_PTR (x[b], x[c]);
     PLL_SWAP_INT (oi[b], oi[c]);
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

