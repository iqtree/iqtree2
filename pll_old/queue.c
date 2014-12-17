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
 * @file queue.c
 */
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
