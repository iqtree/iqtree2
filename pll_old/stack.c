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
 * @file stack.c
 * @brief Generic stack implementation
 *
 * Detailed description to appear soon.
 */
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
  void * item;
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

